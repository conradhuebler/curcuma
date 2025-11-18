/*
 * <Collective Variable Factory>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (November 2025)
 *
 * CV FACTORY - DESIGN PATTERN
 * ============================
 *
 * MOTIVATION:
 * ----------
 * Creating CV objects manually is verbose:
 *
 * ```cpp
 * // Verbose (without factory)
 * std::unique_ptr<CollectiveVariable> cv;
 * if (type == "distance") {
 *     cv = std::make_unique<CV_Distance>(atom_i, atom_j);
 * } else if (type == "angle") {
 *     cv = std::make_unique<CV_Angle>(atom_i, atom_j, atom_k);
 * } else if ...
 * ```
 *
 * With factory:
 * ```cpp
 * auto cv = CVFactory::create("distance", {0, 10});
 * ```
 *
 * DESIGN PATTERN:
 * --------------
 * We use the **Factory Pattern** to decouple CV creation from usage.
 *
 * Benefits:
 * 1. Centralized creation logic
 * 2. Easy to add new CV types (only modify factory)
 * 3. Can create CVs from strings (useful for JSON/command-line input)
 * 4. Type safety (compile-time checks)
 *
 * REFERENCES:
 * ----------
 * [1] Gamma, E. et al. (1994). Design Patterns: Elements of Reusable
 *     Object-Oriented Software. Addison-Wesley. (Factory Pattern, p. 107)
 *
 * [2] PLUMED uses a similar pattern for CV creation (ActionRegister class)
 */

#pragma once

#include "collective_variable.h"
#include "cv_distance.h"
#include "cv_angle.h"
#include "cv_dihedral.h"
#include "cv_gyration.h"
#include "cv_coordination.h"

#include <memory>
#include <string>
#include <stdexcept>
#include <algorithm>

namespace CV {

/**
 * @class CVFactory
 * @brief Factory for creating collective variable objects
 *
 * USAGE:
 * -----
 * ```cpp
 * // Create distance CV
 * auto dist_cv = CVFactory::create("distance", {0, 10});
 *
 * // Create angle CV
 * auto angle_cv = CVFactory::create("angle", {5, 0, 10});
 *
 * // Create dihedral CV
 * auto dih_cv = CVFactory::create("dihedral", {1, 2, 3, 4});
 *
 * // Create radius of gyration (all atoms)
 * auto rg_cv = CVFactory::create("gyration", {});
 *
 * // Create coordination number
 * auto coord_cv = CVFactory::create("coordination", {0}, {1, 5, 9});
 * ```
 *
 * EXTENSIBILITY:
 * -------------
 * To add a new CV type:
 * 1. Create new header cv_mynewcv.h
 * 2. Implement CV_MyNewCV : public CollectiveVariable
 * 3. Add #include "cv_mynewcv.h" to this file
 * 4. Add case to CVFactory::create()
 *
 * @author Claude (Anthropic) & Conrad Hübler
 * @date November 2025
 */
class CVFactory {
public:
    /**
     * @brief Create a collective variable from type string
     *
     * @param type CV type name (case-insensitive):
     *             "distance", "angle", "dihedral", "gyration", "coordination", "rmsd"
     * @param atoms Atom indices (meaning depends on CV type)
     * @param atoms_B Second group of atoms (for coordination CV only)
     * @return Unique pointer to CV object
     *
     * ATOM INDICES CONVENTION:
     * - Distance: {i, j} (2 atoms)
     * - Angle: {i, j, k} (3 atoms, j is vertex)
     * - Dihedral: {i, j, k, l} (4 atoms)
     * - Gyration: {} (empty = all atoms) or {subset of atoms}
     * - Coordination: atoms = group A, atoms_B = group B
     * - RMSD: {} (empty = all atoms) or {subset}
     *
     * ERRORS:
     * - Throws std::invalid_argument if type is unknown
     * - Throws std::invalid_argument if wrong number of atoms
     */
    static std::unique_ptr<CollectiveVariable> create(
        const std::string& type,
        const std::vector<int>& atoms = {},
        const std::vector<int>& atoms_B = {})
    {
        // Convert type to lowercase for case-insensitive matching
        std::string type_lower = type;
        std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        // Distance CV
        if (type_lower == "distance" || type_lower == "dist" || type_lower == "d") {
            if (atoms.size() != 2) {
                throw std::invalid_argument("Distance CV requires exactly 2 atoms, got " +
                                            std::to_string(atoms.size()));
            }
            return std::make_unique<CV_Distance>(atoms[0], atoms[1]);
        }

        // Angle CV
        else if (type_lower == "angle" || type_lower == "ang" || type_lower == "a") {
            if (atoms.size() != 3) {
                throw std::invalid_argument("Angle CV requires exactly 3 atoms, got " +
                                            std::to_string(atoms.size()));
            }
            return std::make_unique<CV_Angle>(atoms[0], atoms[1], atoms[2]);
        }

        // Dihedral CV
        else if (type_lower == "dihedral" || type_lower == "torsion" || type_lower == "dihed") {
            if (atoms.size() != 4) {
                throw std::invalid_argument("Dihedral CV requires exactly 4 atoms, got " +
                                            std::to_string(atoms.size()));
            }
            return std::make_unique<CV_Dihedral>(atoms[0], atoms[1], atoms[2], atoms[3]);
        }

        // Radius of Gyration CV
        else if (type_lower == "gyration" || type_lower == "rg" || type_lower == "rog") {
            auto cv = std::make_unique<CV_Gyration>();
            if (!atoms.empty()) {
                cv->setAtoms(atoms);
            }
            // If atoms is empty, CV will use all atoms
            return cv;
        }

        // Coordination Number CV
        else if (type_lower == "coordination" || type_lower == "coord" || type_lower == "cn") {
            auto cv = std::make_unique<CV_Coordination>();
            if (atoms.empty()) {
                throw std::invalid_argument("Coordination CV requires at least group A atoms");
            }
            cv->setGroupA(atoms);
            if (!atoms_B.empty()) {
                cv->setGroupB(atoms_B);
            }
            // If atoms_B is empty, CV will compute self-coordination
            return cv;
        }

        // RMSD CV (future implementation)
        else if (type_lower == "rmsd") {
            throw std::runtime_error("RMSD CV not yet implemented in factory (use existing RMSD code)");
            // TODO: Implement CV_RMSD wrapper around existing RMSDDriver
        }

        // Unknown type
        else {
            throw std::invalid_argument("Unknown CV type: '" + type +
                                        "'. Valid types: distance, angle, dihedral, gyration, coordination");
        }
    }

    /**
     * @brief Create CV from JSON configuration
     *
     * @param config JSON object with CV specification
     * @return Unique pointer to CV object
     *
     * EXPECTED JSON FORMAT:
     * ```json
     * {
     *   "type": "distance",
     *   "atoms": [0, 10],
     *   "name": "d_Na_O"  // optional
     * }
     * ```
     *
     * or for coordination:
     * ```json
     * {
     *   "type": "coordination",
     *   "atoms": [0],       // group A
     *   "atoms_B": [1,5,9], // group B
     *   "cutoff": 3.5,
     *   "n": 6,
     *   "m": 12
     * }
     * ```
     */
    static std::unique_ptr<CollectiveVariable> createFromJson(const json& config) {
        // Extract type
        if (!config.contains("type")) {
            throw std::invalid_argument("CV config must contain 'type' field");
        }
        std::string type = config["type"];

        // Extract atoms
        std::vector<int> atoms;
        if (config.contains("atoms")) {
            atoms = config["atoms"].get<std::vector<int>>();
        }

        std::vector<int> atoms_B;
        if (config.contains("atoms_B")) {
            atoms_B = config["atoms_B"].get<std::vector<int>>();
        }

        // Create CV
        auto cv = create(type, atoms, atoms_B);

        // Set optional name
        if (config.contains("name")) {
            cv->setName(config["name"]);
        }

        // Type-specific configuration
        if (type == "coordination" || type == "coord") {
            auto* coord_cv = dynamic_cast<CV_Coordination*>(cv.get());
            if (coord_cv) {
                if (config.contains("cutoff")) {
                    coord_cv->setCutoff(config["cutoff"]);
                }
                if (config.contains("n") && config.contains("m")) {
                    coord_cv->setExponents(config["n"], config["m"]);
                }
                if (config.contains("pbc")) {
                    coord_cv->setUsePBC(config["pbc"]);
                }
            }
        }

        return cv;
    }

    /**
     * @brief Get list of all available CV types
     *
     * @return Vector of CV type names
     */
    static std::vector<std::string> getAvailableTypes() {
        return {
            "distance",
            "angle",
            "dihedral",
            "gyration",
            "coordination"
            // "rmsd"  // TODO: Add when CV_RMSD is implemented
        };
    }

    /**
     * @brief Check if a CV type is available
     *
     * @param type CV type name (case-insensitive)
     * @return True if type is valid
     */
    static bool isValidType(const std::string& type) {
        std::string type_lower = type;
        std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        auto types = getAvailableTypes();
        return std::find(types.begin(), types.end(), type_lower) != types.end();
    }

    /**
     * @brief Print help message showing available CV types
     */
    static void printHelp() {
        std::cout << "Available Collective Variable Types:\n";
        std::cout << "====================================\n\n";

        std::cout << "1. DISTANCE (atoms: i, j)\n";
        std::cout << "   - Euclidean distance between two atoms\n";
        std::cout << "   - Usage: CVFactory::create(\"distance\", {0, 10})\n\n";

        std::cout << "2. ANGLE (atoms: i, j, k)\n";
        std::cout << "   - Valence angle at vertex j\n";
        std::cout << "   - Usage: CVFactory::create(\"angle\", {5, 0, 10})\n\n";

        std::cout << "3. DIHEDRAL (atoms: i, j, k, l)\n";
        std::cout << "   - Torsion angle around bond j-k\n";
        std::cout << "   - Usage: CVFactory::create(\"dihedral\", {1, 2, 3, 4})\n\n";

        std::cout << "4. GYRATION (atoms: subset or empty for all)\n";
        std::cout << "   - Radius of gyration (compactness)\n";
        std::cout << "   - Usage: CVFactory::create(\"gyration\", {})\n\n";

        std::cout << "5. COORDINATION (atoms: group A, atoms_B: group B)\n";
        std::cout << "   - Coordination number with smooth switching\n";
        std::cout << "   - Usage: CVFactory::create(\"coordination\", {0}, {1,5,9})\n\n";

        std::cout << "See docs/COLLECTIVE_VARIABLES.md for detailed documentation.\n";
    }
};

} // namespace CV
