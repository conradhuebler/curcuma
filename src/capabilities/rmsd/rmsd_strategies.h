/*
 * <Alignment Strategy Pattern for RMSD calculations>
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
 * Claude Generated - Strategy Pattern implementation for 7 alignment methods
 */

#pragma once

#include "rmsd_assignment.h"
#include "rmsd_costmatrix.h"
#include "src/core/molecule.h"
#include <map>
#include <memory>

// Claude Generated - Alignment method enumeration for better readability
enum class AlignmentMethod : int {
    INCREMENTAL = 1, // Legacy incremental alignment with threading
    TEMPLATE = 2, // Fragment-based template + Kuhn-Munkres
    HEAVY_TEMPLATE = 3, // Heavy atom template (legacy)
    ATOM_TEMPLATE = 4, // Selected atom template + Kuhn-Munkres
    INERTIA = 5, // Inertia-based + Kuhn-Munkres (template-free)
    MOLALIGN = 6, // External MolAlign integration
    DISTANCE_TEMPLATE = 7, // Distance-based template + Kuhn-Munkres
    PREDEFINED_ORDER = 10 // Predefined order (special case)
};

// Forward declarations
class RMSDDriver; // To avoid circular dependency
struct AlignmentConfig; // Forward declaration - defined in rmsd.h

/**
 * @brief Result of alignment calculation
 * Claude Generated - Encapsulates alignment results
 */
struct AlignmentResult {
    double rmsd = 0.0; // Calculated RMSD
    double rmsd_raw = 0.0; // Raw RMSD before optimization
    std::vector<int> reorder_rules; // Final reordering rules
    Molecule reference_aligned; // Aligned reference molecule
    Molecule target_aligned; // Aligned target molecule
    Molecule target_reordered; // Reordered target molecule
    bool success = false; // Whether alignment succeeded
    std::string error_message; // Error message if failed

    AlignmentResult() = default;
};

/**
 * @brief Abstract base class for alignment strategies
 * Claude Generated - Strategy pattern interface
 */
class AlignmentStrategy {
public:
    virtual ~AlignmentStrategy() = default;

    /**
     * @brief Perform molecular alignment using this strategy
     * @param driver Pointer to RMSDDriver for access to molecules and utilities
     * @param config Configuration for this alignment
     * @return Result of alignment calculation
     */
    virtual AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) = 0;

    /**
     * @brief Get human-readable name of this strategy
     */
    virtual std::string getName() const = 0;

    /**
     * @brief Get method enum for this strategy
     */
    virtual AlignmentMethod getMethod() const = 0;
};

/**
 * @brief Method 1: Incremental Alignment (Legacy)
 * Claude Generated - Refactored from ReorderIncremental()
 */
class IncrementalAlignmentStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "Incremental Alignment without Kuhn-Munkres (Legacy)"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::INCREMENTAL; }
};

/**
 * @brief Method 2: Template Approach + Kuhn-Munkres
 * Claude Generated - Refactored from TemplateReorder()
 */
class TemplateAlignmentStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "Prealignment with Template approach and Permutation with Kuhn-Munkres"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::TEMPLATE; }
};

/**
 * @brief Method 3: Heavy Atom Template (Legacy)
 * Claude Generated - Refactored from HeavyTemplate()
 */
class HeavyTemplateStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "Subspace: Prealignment with incremental template approach using all non-hydrogen elements and Permutation with Kuhn-Munkres (Legacy)"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::HEAVY_TEMPLATE; }
};

/**
 * @brief Method 4: Selected Atom Template + Kuhn-Munkres
 * Claude Generated - Refactored from AtomTemplate()
 */
class AtomTemplateStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "Subspace: Prealignment with incremental template approach using selected elements and Permutation with Kuhn-Munkres"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::ATOM_TEMPLATE; }
};

/**
 * @brief Method 5: Inertia-based Alignment + Kuhn-Munkres
 * Claude Generated - Refactored from TemplateFree()
 */
class InertiaAlignmentStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "Inertia: Prealignment with according to Moment of Inertia and Permutation with Kuhn-Munkres"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::INERTIA; }
};

/**
 * @brief Method 6: External MolAlign Integration
 * Claude Generated - Refactored from MolAlignLib()
 */
class MolAlignStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "MolAlign: Using external molalign for permutation"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::MOLALIGN; }
};

/**
 * @brief Method 7: Distance Template + Kuhn-Munkres
 * Claude Generated - Refactored from DistanceTemplate()
 */
class DistanceTemplateStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "Distance template: Prealignment with incremental template approach using distance criterion and Permutation with Kuhn-Munkres"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::DISTANCE_TEMPLATE; }
};

/**
 * @brief Method 10: Predefined Order (Special case)
 * Claude Generated - Handles predefined reorder rules
 */
class PredefinedOrderStrategy : public AlignmentStrategy {
public:
    AlignmentResult align(RMSDDriver* driver, const AlignmentConfig& config) override;
    std::string getName() const override { return "Predefined Order"; }
    AlignmentMethod getMethod() const override { return AlignmentMethod::PREDEFINED_ORDER; }
};

/**
 * @brief Factory for creating alignment strategies
 * Claude Generated - Strategy factory pattern
 */
class AlignmentStrategyFactory {
public:
    /**
     * @brief Create strategy based on alignment method
     * @param method Alignment method enum
     * @return Unique pointer to strategy, nullptr if invalid method
     */
    static std::unique_ptr<AlignmentStrategy> createStrategy(AlignmentMethod method);

    /**
     * @brief Create strategy based on method ID (legacy compatibility)
     * @param method_id Method identifier (1-7, 10)
     * @return Unique pointer to strategy, nullptr if invalid method
     */
    static std::unique_ptr<AlignmentStrategy> createStrategy(int method_id);

    /**
     * @brief Get all available strategies with their method IDs and names
     */
    static std::map<int, std::string> getAvailableStrategies();
};