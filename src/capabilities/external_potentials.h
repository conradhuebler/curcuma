/*
 * External potentials for SimpleMD — a configured, persistent bias.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - SimpleMD::applyExternalForces() takes a per-atom force
 * matrix and clears it after one step. That is an injection: the caller has to
 * re-apply it every step, nothing about it is written down, and a run cannot be
 * reproduced from its configuration. These are potentials instead. They live in
 * the controller, they are evaluated inside the force loop the way the walls are,
 * and they can be replaced while the run is going.
 *
 * Atom sets are named with the existing selection grammar (FragString2Indicies),
 * so "F2" means the second fragment here exactly as everywhere else.
 *
 * Units follow SimpleMD's own: geometry in Angstrom, gradient in Eh/Angstrom,
 * energy in Eh. A force constant is therefore Eh/Angstrom^2 and a constant force
 * Eh/Angstrom.
 */

#pragma once

#include "src/core/global.h"
#include "src/core/molecule.h"

#include <array>
#include <string>
#include <vector>

#include "json.hpp"
using json = nlohmann::json;

namespace curcuma {

/**
 * @brief One configured external potential.
 *
 * @a work accumulates the work the potential has done on the system,
 * sum over steps of F . dr. That is what turns a steered pull from something
 * visible into something quantitative, and it is the number a Jarzynski or Crooks
 * estimate is built from.
 */
struct ExternalPotential {
    enum class Kind {
        ConstantForce,     ///< a fixed force on every atom of the set
        CentroidHarmonic,  ///< the set's centroid restrained to a point
        DistanceHarmonic   ///< the two sets' centroids restrained to a distance
    };

    Kind kind = Kind::ConstantForce;
    std::string label;                     ///< for the log and the results
    std::string selection;                 ///< the expression as written, kept for the record
    std::string selection_b;
    std::vector<int> atoms;                ///< resolved, zero-based
    std::vector<int> atoms_b;

    std::array<double, 3> direction { { 0.0, 0.0, 0.0 } };  ///< normalised on parse
    double magnitude = 0.0;                                 ///< Eh/Angstrom
    std::array<double, 3> target { { 0.0, 0.0, 0.0 } };     ///< Angstrom
    double k = 0.0;                                         ///< Eh/Angstrom^2
    double r0 = 0.0;                                        ///< Angstrom

    double energy = 0.0;   ///< of the last evaluation, Eh
    double value = 0.0;    ///< the coordinate it acts on: distance, or |R_c - T|
    double work = 0.0;     ///< accumulated F . dr since the run started, Eh
};

/**
 * @brief Read a list of potentials from the controller.
 *
 * @p list is the JSON array the `external_potentials` parameter holds. Selection
 * strings are resolved against @p molecule, which is also where an unknown or
 * empty selection is caught -- a bias on nothing is a configuration error, not a
 * no-op to be discovered later. Returns an empty vector and sets @p error on the
 * first entry that does not parse.
 */
std::vector<ExternalPotential> parseExternalPotentials(const json& list,
    const Molecule& molecule, std::string* error);

/**
 * @brief Add the potentials' energy and gradient at @p geometry.
 *
 * The gradient (dE/dr, Eh/Angstrom) is accumulated into @p gradient, matching the
 * sign convention of the wall potentials. When @p previous has the same shape as
 * @p geometry, each potential's accumulated work is advanced by F . dr over the
 * step just taken. Returns the summed energy in Eh.
 */
double applyExternalPotentials(std::vector<ExternalPotential>& potentials,
    const Geometry& geometry, const Geometry& previous, Geometry& gradient);

/** @brief The potentials, their settings and what they have done, for Results(). */
json describeExternalPotentials(const std::vector<ExternalPotential>& potentials);

} // namespace curcuma
