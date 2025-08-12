/*
 * EnergyCalculator Refactoring Helper Functions
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated: Helper functions to simplify repetitive switch statements
 */

#pragma once

#include "energycalculator_enums.h"
#include <functional>

/**
 * @brief Action types for different EnergyCalculator operations
 */
enum class ActionType {
    INITIALIZE, ///< Initialize interfaces and allocate resources
    CLEANUP, ///< Clean up and delete allocated resources
    SET_MOLECULE ///< Set molecule data in interfaces
};

/**
 * @brief Unified method dispatcher for EnergyCalculator operations
 *
 * This class replaces the repetitive switch statements in constructor,
 * destructor, and setMolecule with a single unified dispatch system.
 */
class MethodDispatcher {
public:
    /**
     * @brief Execute action based on method type
     * @param calc Reference to EnergyCalculator instance
     * @param method_type Method type to dispatch to
     * @param action Type of action to perform
     * @param mol Optional molecule parameter for SET_MOLECULE action
     */
    static void dispatch(EnergyCalculator* calc, MethodType method_type,
        ActionType action, const Mol* mol = nullptr);

private:
    // Action handlers for each method type
    static void handleForceField(EnergyCalculator* calc, ActionType action, const Mol* mol);
    static void handleTBLite(EnergyCalculator* calc, ActionType action, const Mol* mol);
    static void handleXTB(EnergyCalculator* calc, ActionType action, const Mol* mol);
    static void handleUlysses(EnergyCalculator* calc, ActionType action, const Mol* mol);
    static void handleDFTD3(EnergyCalculator* calc, ActionType action, const Mol* mol);
    static void handleDFTD4(EnergyCalculator* calc, ActionType action, const Mol* mol);
    static void handleEHT(EnergyCalculator* calc, ActionType action, const Mol* mol);

    // Helper functions
    static bool checkCompilationFlag(const std::string& flag);
    static void exitWithMessage(const std::string& message);
};