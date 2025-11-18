/*
 * <Collective Variable Framework for Enhanced Sampling>
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
 * ---
 * Claude Generated (November 2025)
 *
 * COLLECTIVE VARIABLES FOR METADYNAMICS
 *
 * This module implements a polymorphic framework for collective variables (CVs)
 * used in metadynamics simulations. Inspired by PLUMED's architecture [1].
 *
 * THEORETICAL BACKGROUND:
 * ----------------------
 * In metadynamics, a bias potential V(s,t) is constructed as a function of
 * collective variables s(r), where r are atomic Cartesian coordinates:
 *
 *     V(s,t) = Σ_i w_i * exp(-(s - s_i)² / (2σ²))
 *
 * The bias force on atom i is computed via the chain rule:
 *
 *     F_i^bias = -(∂V/∂s) * (∂s/∂r_i)
 *
 * where:
 *   - ∂V/∂s: Derivative of bias potential w.r.t. CV (from MTD algorithm)
 *   - ∂s/∂r_i: Gradient of CV w.r.t. atomic positions (Jacobian, from CV class)
 *
 * KEY REQUIREMENTS FOR CVs:
 * ------------------------
 * 1. Differentiability: ∇s(r) must be well-defined
 * 2. Physical relevance: Captures slow degrees of freedom
 * 3. Low-dimensionality: Typically 1-3 CVs (curse of dimensionality)
 * 4. Computational efficiency: Fast evaluation during MD
 *
 * REFERENCES:
 * ----------
 * [1] Bonomi, M. et al. (2009). PLUMED: A portable plugin for free-energy
 *     calculations with molecular dynamics. Comp. Phys. Comm. 180, 1961-1972.
 *     DOI: 10.1016/j.cpc.2009.05.011
 *
 * [2] Laio, A. & Parrinello, M. (2002). Escaping free-energy minima.
 *     Proc. Natl. Acad. Sci. USA 99, 12562-12566.
 *     DOI: 10.1073/pnas.202427399
 *
 * [3] Tribello, G. A. et al. (2014). PLUMED 2: New feathers for an old bird.
 *     Comp. Phys. Comm. 185, 604-613.
 *     DOI: 10.1016/j.cpc.2013.09.018
 *
 * DOCUMENTATION:
 * -------------
 * See docs/COLLECTIVE_VARIABLES.md for comprehensive documentation including:
 *   - Mathematical formulations with gradients
 *   - Use cases and examples
 *   - Numerical stability considerations
 *   - Performance optimization strategies
 */

#pragma once

#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/tools/geometry.h"

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace CV {

/**
 * @enum CVType
 * @brief Enumeration of all implemented collective variable types
 *
 * Each CV type corresponds to a specific derived class implementing
 * the CollectiveVariable interface.
 */
enum class CVType {
    RMSD,           ///< Root-Mean-Square Deviation from reference structure
    Distance,       ///< Distance between two atoms or centers of mass
    Angle,          ///< Angle formed by three atoms (i-j-k)
    Dihedral,       ///< Dihedral/torsion angle formed by four atoms (i-j-k-l)
    Gyration,       ///< Radius of gyration (compactness measure)
    Coordination,   ///< Coordination number (smooth switching function)
    PathCV,         ///< Path collective variable (s,z) [future]
    Custom          ///< User-defined CV via lambda function [future]
};

/**
 * @brief Convert CVType enum to human-readable string
 */
inline std::string CVTypeToString(CVType type) {
    switch (type) {
    case CVType::RMSD: return "RMSD";
    case CVType::Distance: return "Distance";
    case CVType::Angle: return "Angle";
    case CVType::Dihedral: return "Dihedral";
    case CVType::Gyration: return "RadiusOfGyration";
    case CVType::Coordination: return "Coordination";
    case CVType::PathCV: return "PathCV";
    case CVType::Custom: return "Custom";
    default: return "Unknown";
    }
}

/**
 * @class CollectiveVariable
 * @brief Abstract base class for all collective variables
 *
 * DESIGN PHILOSOPHY:
 * -----------------
 * This class defines the interface that all CVs must implement. The polymorphic
 * design allows for easy extension with new CV types without modifying existing
 * metadynamics code (Open/Closed Principle).
 *
 * KEY METHODS:
 * -----------
 * - calculate(): Compute current CV value s(r) from molecular geometry
 * - gradient(): Compute Jacobian ∂s/∂r_i (N_atoms × 3 matrix)
 * - type(): Return CV type identifier
 *
 * USAGE PATTERN:
 * -------------
 * ```cpp
 * std::unique_ptr<CollectiveVariable> cv = CVFactory::create("distance", {0, 10});
 * double s = cv->calculate(molecule);
 * Geometry grad = cv->gradient(molecule);  // For bias force calculation
 * ```
 *
 * IMPORTANT NOTES:
 * ---------------
 * - All gradients are returned in Cartesian coordinates (Å or degrees)
 * - Gradient matrix has shape (N_atoms × 3), with zeros for uninvolved atoms
 * - Thread-safety: CV objects are not thread-safe; create one per thread if needed
 *
 * @author Claude (Anthropic) & Conrad Hübler
 * @date November 2025
 */
class CollectiveVariable {
public:
    virtual ~CollectiveVariable() = default;

    /**
     * @brief Calculate the current value of the collective variable
     *
     * @param mol Molecule object containing current atomic positions
     * @return Current CV value (units depend on CV type: Å, degrees, etc.)
     *
     * PERFORMANCE:
     * This is called every MD step, so implementations must be efficient.
     * Typical cost: O(N_atoms) for most CVs, O(N_atoms²) for RMSD.
     *
     * THREAD-SAFETY:
     * Not thread-safe. Use separate CV instances per thread.
     */
    virtual double calculate(const Molecule& mol) = 0;

    /**
     * @brief Calculate the gradient ∂s/∂r_i (Jacobian matrix)
     *
     * @param mol Molecule object containing current atomic positions
     * @return Geometry matrix (N_atoms × 3) with gradient components
     *
     * MATHEMATICAL BACKGROUND:
     * The gradient is used to compute bias forces via chain rule:
     *     F_i^bias = -(∂V/∂s) * (∂s/∂r_i)
     *
     * NUMERICAL STABILITY:
     * - Implementations must handle singularities (e.g., d→0 for distances)
     * - Add small epsilon (10^-8) to denominators where needed
     * - Use atan2() for angle calculations to handle all quadrants
     *
     * TESTING:
     * Gradient correctness is verified via finite differences:
     *     ∂s/∂x_i ≈ (s(x_i + ε) - s(x_i - ε)) / (2ε)
     *
     * RETURN VALUE:
     * - Matrix has zeros for atoms not involved in the CV
     * - Units: Same as CV value divided by Å (e.g., unitless for RMSD, Å^-1 for distance)
     */
    virtual Geometry gradient(const Molecule& mol) = 0;

    /**
     * @brief Get the type of this collective variable
     */
    virtual CVType type() const = 0;

    /**
     * @brief Get human-readable description of this CV
     *
     * Example: "Distance(0,10)" or "Dihedral(φ, atoms 5-6-7-8)"
     */
    virtual std::string description() const = 0;

    /**
     * @brief Get name/label for this CV (user-defined)
     */
    inline std::string name() const { return m_name; }

    /**
     * @brief Set name/label for this CV
     */
    inline void setName(const std::string& name) { m_name = name; }

    /**
     * @brief Get list of atoms involved in this CV
     *
     * @return Vector of 0-based atom indices
     *
     * NOTE: Empty vector means "all atoms" (e.g., for radius of gyration)
     */
    inline const std::vector<int>& getAtoms() const { return m_atoms; }

    /**
     * @brief Set atoms involved in this CV
     *
     * @param atoms Vector of 0-based atom indices
     *
     * VALIDATION:
     * Derived classes should validate that the number of atoms is correct
     * (e.g., Distance requires 2 atoms, Angle requires 3, etc.)
     */
    virtual void setAtoms(const std::vector<int>& atoms) {
        m_atoms = atoms;
    }

    /**
     * @brief Check if CV is periodic (e.g., dihedral angles wrap at ±180°)
     *
     * This is important for:
     * - Correctly computing differences between CV values
     * - FES reconstruction (histogram binning)
     * - Grid-based acceleration
     */
    virtual bool isPeriodic() const { return false; }

    /**
     * @brief Get period for periodic CVs (e.g., 360° for angles)
     *
     * @return Period value, or 0.0 for non-periodic CVs
     */
    virtual double getPeriod() const { return 0.0; }

    /**
     * @brief Compute difference between two CV values, accounting for periodicity
     *
     * For non-periodic CVs: simply (a - b)
     * For periodic CVs: smallest distance on periodic domain
     *
     * Example (dihedral):
     *   diff(+170°, -170°) = -20° (not +340°)
     */
    virtual double periodicDifference(double a, double b) const {
        if (!isPeriodic()) {
            return a - b;
        }
        double period = getPeriod();
        double diff = a - b;
        // Wrap to [-period/2, +period/2]
        while (diff > period / 2.0) diff -= period;
        while (diff < -period / 2.0) diff += period;
        return diff;
    }

protected:
    std::vector<int> m_atoms;  ///< Atom indices involved in this CV (0-based)
    std::string m_name;         ///< User-defined label for this CV
};

/**
 * @brief Smart pointer type for collective variables
 *
 * MEMORY MANAGEMENT:
 * Use unique_ptr for ownership, shared_ptr if CVs need to be shared across threads.
 */
using CVPtr = std::unique_ptr<CollectiveVariable>;

} // namespace CV
